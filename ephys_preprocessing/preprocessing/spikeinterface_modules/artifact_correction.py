from __future__ import annotations

import numpy as np
from typing import Union, Optional

from spikeinterface.core.core_tools import define_function_from_class
from spikeinterface.core import BaseRecording, get_chunk_with_margin
from spikeinterface.preprocessing.basepreprocessor import BasePreprocessor, BasePreprocessorSegment


class ArtifactCorrectionRecording(BasePreprocessor):
    """
    Recording preprocessor that corrects artifacts by replacing artifact timepoints
    with the mean of preceding data samples. This implementation exactly matches
    the behavior of the memory-mapped version.
    
    Parameters
    ----------
    recording : BaseRecording
        The recording to be corrected
    artifact_times : numpy.ndarray
        Array of artifact timestamps in seconds
    window_ms : float, default: 3.0
        Time window in milliseconds for artifact correction. This window size is used
        both for the replacement window and for computing the mean of preceding data.
    dtype : dtype, default: None
        The dtype of the returned traces. If None, the dtype is maintained.
        
    Returns
    -------
    corrected_recording : ArtifactCorrectionRecording
        The artifact-corrected recording
    """
    
    def __init__(self, 
                 recording: BaseRecording, 
                 artifact_times: np.ndarray, 
                 window_ms: float = 3.0, 
                 dtype: Optional[np.dtype] = None):
        BasePreprocessor.__init__(self, recording, dtype=dtype)
        
        artifact_times = np.asarray(artifact_times)
        
        # Convert artifact times to sample indices
        fs = self.get_sampling_frequency()
        indices = np.round(artifact_times * fs).astype(int)
        window_samples = int(window_ms * fs / 1000)
        
        # Clip indices to valid range (same as memmap version)
        indices = np.clip(indices, 0, recording.get_num_samples() - 1)
        
        # Create chunks of indices for the correction window
        correction_indices = np.array([np.arange(i, i + window_samples) for i in indices])
        
        for parent_segment in recording._recording_segments:
            rec_segment = ArtifactCorrectionRecordingSegment(
                parent_segment,
                correction_indices,
                window_samples,
                dtype
            )
            self.add_recording_segment(rec_segment)

        self._kwargs = dict(
            recording=recording,
            artifact_times=artifact_times,
            window_ms=window_ms,
            dtype=dtype
        )


class ArtifactCorrectionRecordingSegment(BasePreprocessorSegment):
    def __init__(self, parent_recording_segment, correction_indices, window_samples, dtype):
        BasePreprocessorSegment.__init__(self, parent_recording_segment)
        self.correction_indices = correction_indices
        self.window_samples = window_samples
        self.dtype = dtype

    def get_traces(self, start_frame, end_frame, channel_indices):
        if start_frame is None:
            start_frame = 0
        if end_frame is None:
            end_frame = self.get_num_samples()
        if channel_indices is None:
            channel_indices = slice(None)

        # Get traces with extra margin for computing means
        traces, left_margin, right_margin = get_chunk_with_margin(
            self.parent_recording_segment,
            start_frame,
            end_frame,
            channel_indices,
            self.window_samples,  # Use window_samples as margin
            add_zeros=True  # Add zeros padding if needed
        )
        
        # Find artifacts that fall within the requested time window
        relevant_artifacts = []
        for artifact_indices in self.correction_indices:
            # Check if any part of the artifact window overlaps with the requested chunk
            if (artifact_indices[0] < end_frame) and (artifact_indices[-1] >= start_frame):
                relevant_artifacts.append(artifact_indices)
        
        if relevant_artifacts:
            # Make a copy since we'll modify the data
            traces = traces.copy()
            
            for artifact_indices in relevant_artifacts:
                # Calculate indices relative to the chunk
                chunk_relative_indices = artifact_indices - start_frame + left_margin
                
                # Only process indices that fall within the chunk
                valid_mask = (chunk_relative_indices >= 0) & (chunk_relative_indices < traces.shape[0])
                chunk_relative_indices = chunk_relative_indices[valid_mask]
                
                if len(chunk_relative_indices) > 0:
                    # Get mean of data exactly window_samples before the artifact (matching memmap version)
                    before_indices = chunk_relative_indices[0] - self.window_samples
                    if before_indices >= 0:
                        # Calculate mean over the same window size
                        before_window = np.arange(before_indices, before_indices + self.window_samples)
                        ch_means = np.mean(traces[before_window], axis=0)
                        
                        # Replace artifact timepoints with repeated mean values (matching memmap version)
                        for i in range(len(chunk_relative_indices)):
                            traces[chunk_relative_indices[i]] = ch_means

        # Remove the margins we added
        if right_margin > 0:
            traces = traces[left_margin:-right_margin]
        else:
            traces = traces[left_margin:]

        return traces.astype(self.dtype) if self.dtype else traces


# Function for API 
artifact_correction = define_function_from_class(source_class=ArtifactCorrectionRecording, name="artifact_correction")