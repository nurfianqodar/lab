import requests
import argparse
from pydantic import BaseModel, Field
from typing import Optional
import sys

class Config(BaseModel):
    protein_id: str
    format_file: Optional[str] = Field(default="pdb")

    

def config() -> Config:
    parser = argparse.ArgumentParser("Download protein structure file")
    parser.add_argument("--id", "-i", type=str, help="Protein id. ex: 1AKI")
    parser.add_argument("--format", "-f", type=str, help="File Format [pdb]. ex: pdb", default="pdb")
    args = parser.parse_args()
    return Config(protein_id=args.id, format_file=args.format)


def downloader(id: str, format: str):
    file_name = f"{id}.{format}" 
    print(f"Downloading {file_name} file from PDB")
    try:
        url = f"https://files.rcsb.org/download/{file_name}"
        response = requests.get(url)
        code = response.status_code
        if code != 200:
            raise ProcessLookupError
        else:
            content  = response.content
            with open(file_name , "wb") as file:
                file.write(content)
            print(f"Success download {file_name}")
    except Exception as e:
        print(e)
        sys.exit(1)


conf = config()
downloader(conf.protein_id, conf.format_file)




