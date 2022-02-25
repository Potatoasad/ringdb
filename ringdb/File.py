import os
import subprocess
import h5py

class File:
    def __init__(self, rel_path):
        self.path = rel_path
        
    @classmethod
    def from_url(cls, url, save_folder, new_filename=None):
        # Remove trailing "/" from path
        if save_folder[-1] == "/":
            save_folder = save_folder[0:-1]
        
        # Check if folder one step up exists:
        folder_up = save_folder.split('/')[-2]
        if folder_up != ".":
            if not os.path.exists(folder_up):
                print(f"making folder {folder_up} since it doesn't exist")
                subprocess.run(["mkdir", folder_up])
                
        thefilepath = f"{save_folder}/{url.split('/')[-1]}"
        # Downloading the file we have into the folder
        print(f"Downloading file from {url}")
        subprocess.run(["wget",url,"-P",save_folder])
        file = cls(thefilepath)
        if new_filename is not None:
            file.rename(new_filename)
            
        return file
    
    def extract_here(self):
        file_type = self.path.split('.')[-1]
        folder = "/".join(self.path.split('/')[0:-1])
        if file_type == 'tar':
            subprocess.run(["tar","-xvf",self.path,"-C",folder])
        if file_type == 'zip':
            subprocess.run(["unzip",self.path,"-d",folder])
    
    @property
    def exists(self):
        return os.path.exists(self.path)
        
    def delete(self):
        if self.path not in ['.', './', '/']:
            subprocess.run(["rm", self.path])
        else:
            print(f"Was about to delete {self.path} , aborted it")
        
    def rename(self,newname):
        filename = self.path.split('/')[-1]
        newpath = '/'.join(self.path.split('/')[0:-1] + [newname])
        subprocess.run(["mv",self.path, newpath])
        self.path = newpath

