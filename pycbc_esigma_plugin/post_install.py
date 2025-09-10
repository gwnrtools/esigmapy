import os
import sys

def main():
    lal_data_path = os.environ.get('LAL_DATA_PATH')
    if not lal_data_path:
        lal_data_path = input('Please enter the path for LAL_DATA_PATH: ')
    
    if lal_data_path:
        try:
            import gwnr
            
            # Path to esigma_utils.py using sys.modules
            esigma_utils_path = os.path.join(os.path.dirname(sys.modules['gwnr'].__file__), 'waveform', 'esigma_utils.py')
            
            # Read the content of esigma_utils.py
            with open(esigma_utils_path, 'r') as file:
                lines = file.readlines()
            
            # Define the lines to insert
            import_os_line = 'import os\n'
            lal_data_path_line = f'os.environ["LAL_DATA_PATH"] = "{lal_data_path}"\n'
            
            # Check if the lines already exist
            import_os_exists = any(line.strip() == import_os_line.strip() for line in lines)
            lal_data_path_exists = any(line.strip() == lal_data_path_line.strip() for line in lines)
            
            if import_os_exists:
                print(f'The line "{import_os_line.strip()}" already exists in esigma_utils.py')
            if lal_data_path_exists:
                print(f'The line "{lal_data_path_line.strip()}" already exists in esigma_utils.py')
            
            if not import_os_exists or not lal_data_path_exists:
                # Find the index where to insert the LAL_DATA_PATH line
                insert_index = None
                for index, line in enumerate(lines):
                    if line.startswith('import'):
                        insert_index = index + 1
                
                if insert_index is not None:
                    # Insert the lines if they don't already exist
                    if not import_os_exists:
                        lines.insert(insert_index, import_os_line)
                    if not lal_data_path_exists:
                        lines.insert(insert_index + 1, lal_data_path_line)
                    
                    # Write the updated content back to esigma_utils.py
                    with open(esigma_utils_path, 'w') as file:
                        file.writelines(lines)
                    
                    print(f'LAL_DATA_PATH has been set to: {lal_data_path} in esigma_utils.py')
                else:
                    print('Could not find a suitable place to insert the LAL_DATA_PATH line.')
        except ImportError:
            print('The gwnr module could not be imported. Ensure it is installed.')
        except Exception as e:
            print(f'An error occurred while modifying esigma_utils.py: {e}')
    else:
        print('LAL_DATA_PATH is not set. Exiting.')

if __name__ == "__main__":
    main()
