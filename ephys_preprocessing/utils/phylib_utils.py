from phylib.io.model import TemplateModel, get_template_params
import numpy as np
import pandas as pd

class ExtendedTemplateModel(TemplateModel):
    def get_template_amplitude(self, cluster_id):
        """Return the maximum amplitude of a template's waveforms across all channels."""
        waveforms = self.get_template(cluster_id).template
        assert waveforms.ndim == 2  # shape: (n_samples, n_channels)
        return (waveforms.max(axis=0) - waveforms.min(axis=0)).max()

    def get_best_channels(self, cluster_id):
        """Return the best channels of a given cluster."""
        template = self.get_template(cluster_id)
        if not template:  # pragma: no cover
            return [0]
        return template.channel_ids

    def get_best_channel(self, cluster_id):
        """Return the best channel id of a given cluster."""
        channel_ids = self.get_best_channels(cluster_id)
        assert channel_ids is not None and len(channel_ids)
        return channel_ids[0]

    def get_channel_shank(self, cluster_id):
        """Return the shank of a cluster's best channel."""
        best_channel_id = self.get_best_channel(cluster_id)
        return self.channel_shanks[best_channel_id]

    def get_probe_depth(self, cluster_id):
        """Return the depth of a cluster."""
        channel_id = self.get_best_channel(cluster_id)
        return 0 if channel_id is None else self.channel_positions[channel_id, 1]

    def get_n_spikes(self, cluster_id):
        """Number of spikes in a given cluster."""
        return self.spike_clusters[self.spike_clusters == cluster_id].shape[0]

    def get_mean_firing_rate(self, cluster_id):
        """Return the mean firing rate of a cluster."""
        n_spikes = self.get_n_spikes(cluster_id)
        return n_spikes / max(1, self.duration)

    def create_metrics_dataframe(self):
        """Create a DataFrame with all metrics for each cluster."""
        cluster_ids = self.cluster_ids
        
        df = pd.DataFrame.from_dict(self.metadata)
        
        # Remove all nan-only rows (ContamPct==100) that come from kilosort removing some clusters at the end
        cols_to_keep = [c for c in df.columns if c not in ['Amplitude', 'ContamPct', 'KSLabel']]
        df_to_drop = df[cols_to_keep]
        nan_row_indices = df_to_drop.index[df_to_drop.isnull().all(1)].tolist()
        df = df.drop(nan_row_indices)

        df['cluster_id'] = cluster_ids
        
        # Calculate all metrics
        df['amp'] = df['cluster_id'].apply(self.get_template_amplitude)
        df['ch'] = df['cluster_id'].apply(self.get_best_channel)
        df['sh'] = df['cluster_id'].apply(self.get_channel_shank)
        df['depth'] = df['cluster_id'].apply(self.get_probe_depth)
        df['n_spikes'] = df['cluster_id'].apply(self.get_n_spikes)
        df['fr'] = df['cluster_id'].apply(self.get_mean_firing_rate)
        
        # TODO: group is always nan in my data but maybe should add the actual functions that leads to this
        df['group'] = np.nan

        # Reorder columns
        df = df.reindex(sorted(df.columns), axis=1)
        cluster_ids_column = df.pop('cluster_id')
        df.insert(0, 'cluster_id', cluster_ids_column)

        return df

    def save_metrics_tsv(self, output_path):
        """Save metrics to a TSV file."""
        df = self.create_metrics_dataframe()
        df.to_csv(output_path, sep='\t', index=False)
        return output_path

def load_phy_model(params_path):
    """Return an ExtendedTemplateModel instance from a path to a params.py file."""
    return ExtendedTemplateModel(**get_template_params(params_path))


if __name__ == '__main__':
    # Load the extended model
    model = load_phy_model('path/to/params.py')

    # Create the dataframe
    df = model.create_metrics_dataframe()

    # Save to TSV
    model.save_metrics_tsv('path/to/cluster_info.tsv')