import json
import os

def generate_json(input_string):
    parts = input_string.split('_')
    Nx, Ny, Re, conv, diff, under_relax_u, under_relax_v = map(str, parts[:7])

    data = {
        "in": {
            "uLid": 1.0,
            "geometry": {
                "x": 1.0,
                "y": 1.0
            },
            "fluid": {
                "Re": int(Re)
            }
        },
        "engine": {
            "mesh": {
                "nodes": {
                    "Nx": int(Nx),
                    "Ny": int(Ny)
                }
            },
            "method": {
                "type": "SIMPLE",
                "tolerance": 1e-5,
                "maxIter": 5000,
                "args": {
                    "under_relaxation": {
                        "u": float(under_relax_u) / (10**(len(under_relax_u)-1)),
                        "v": float(under_relax_v) / (10**(len(under_relax_v)-1)),
                        "p": 0.3
                    },
                    "number_of_sweeps": {
                        "u": 4,
                        "v": 4,
                        "p": 8
                    }
                }
            },
            "schemes": {
                "convection": conv,
                "diffusion": diff
            }
        }
    }

    return data

def save_json(data, filename):
    with open(filename, 'w') as outfile:
        json.dump(data, outfile, indent=4)

if __name__ == "__main__":
    file_names = [
        "40_40_100_UDS_SECOND_07_07",
        "40_40_400_UDS_SECOND_07_07",
        "40_40_1000_UDS_SECOND_07_07",
        "40_40_1000_CDS_SECOND_07_07",
        "40_40_1000_QUICK_SECOND_07_07",
        "80_80_1000_UDS_SECOND_07_07",
        "80_80_1000_CDS_SECOND_07_07",
        "80_80_1000_QUICK_SECOND_07_07",
        "129_129_1000_UDS_SECOND_07_07",
        "129_129_1000_CDS_SECOND_07_07",
        "129_129_1000_QUICK_SECOND_07_07",
        "129_129_1000_UDS_FOURTH_07_07",
        "129_129_1000_CDS_FOURTH_07_07",
        "129_129_1000_QUICK_FOURTH_07_07",
    ]

    for filename in os.listdir('.'):
        if filename.endswith('.json') and filename != 'input.json':
            os.remove(filename)
            print(f"Deleted existing file: {filename}")

    for i, file_name in enumerate(file_names, start=0):
        generated_json = generate_json(file_name)
        output_filename = f"{str(i).zfill(2)}_{file_name}.json"
        save_json(generated_json, output_filename)
        print(f"Generated JSON saved to {output_filename}")
