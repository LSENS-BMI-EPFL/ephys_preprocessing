import os
import re

def update_params_files(root_dir):
    """
    Recursively find and update params.py files, replacing dat_path with r'./temp_wh.dat'
    
    Args:
        root_dir (str): Root directory to start searching from
    """
    # Counter for modified files
    modified_count = 0
    
    # Walk through all directories
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if ('AB133' not in dirpath) or ('AB086' not in dirpath) or ('AB121' not in dirpath):
            if 'kilosort4' not in dirpath:
                continue
        if 'params.py' in filenames:
            file_path = os.path.join(dirpath, 'params.py')
            print(f"Found params.py at: {file_path}")
            
            # Read the file content
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Use regex to find and replace the dat_path line
            # This will match dat_path = 'any/path/here' or dat_path = "any/path/here"
            new_content = re.sub(
                r"dat_path\s*=\s*['\"].*?['\"]",
                "dat_path = r'./temp_wh.dat'",
                content
            )
            
            # Only write if the content has changed
            if new_content != content:
                print(f"Modifying: {file_path}")
                # Create backup
                backup_path = file_path + '.bak'
                os.rename(file_path, backup_path)
                
                # Write new content
                with open(file_path, 'w') as f:
                    f.write(new_content)
                
                modified_count += 1
            else:
                print(f"No changes needed in: {file_path}")
    
    return modified_count

if __name__ == "__main__":
    # Get the current directory or specify your root directory
    root_dir = "/mnt/lsens/analysis/Jules_Lebert/data"
    root_dir = "M:/analysis//Axel_Bisi/data"

    print(f"Starting search in: {root_dir}")
    modified = update_params_files(root_dir)
    print(f"\nComplete! Modified {modified} files.")
    print("Backup files were created with .bak extension")