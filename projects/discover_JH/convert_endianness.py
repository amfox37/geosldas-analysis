import os
import struct

def convert_files_in_directory(input_directory, output_directory, data_type):
    for filename in os.listdir(input_directory):
        input_file = os.path.join(input_directory, filename)
        output_file = os.path.join(output_directory, filename)

        with open(input_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                while True:
                    # Read the data in chunks based on the data_type size
                    data = f_in.read(struct.calcsize(data_type))

                    # Break the loop when there is no more data to read
                    if not data:
                        break

                    # Check if we have enough data to unpack
                    if len(data) != struct.calcsize(data_type):
                        print(f"Skipping file {filename} as it does not contain complete data.")
                        break

                    # Unpack the data based on the input file's endianness
                    value = struct.unpack('>' + data_type, data)[0]  # For big-endian use '>'.

                    # Pack the data in little-endian format and write it to the output file
                    f_out.write(struct.pack('<' + data_type, value))

if __name__ == "__main__":
    input_directory = "../test_data/Metop_B/Y2015/M08"
    output_directory = "../test_data/Metop_B/Y2015/M08_le"
    data_type = 'd'  # 'd' represents a 64-bit double-precision floating-point number

    convert_files_in_directory(input_directory, output_directory, data_type)
    print("Conversion complete.")
