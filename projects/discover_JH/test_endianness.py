import struct

# Define the data type and test values
data_type = 'd'  # 'd' represents a 64-bit double-precision floating-point number
test_value1 = 2015
test_value2 = 1.0

# Read the binary file and extract the data
with open('test_file.bin', 'rb') as file:
    data = file.read()

# Try interpreting the data as little-endian and big-endian doubles
little_endian_value = struct.unpack('<' + data_type, data[:8])[0]
big_endian_value = struct.unpack('>' + data_type, data[:8])[0]

print(little_endian_value, big_endian_value)

# Check which test value matches one of the interpreted values
if test_value1 in [little_endian_value, big_endian_value]:
    print("The binary file contains little-endian 64-bit doubles.")
elif test_value2 in [little_endian_value, big_endian_value]:
    print("The binary file contains big-endian 64-bit doubles.")
else:
    print("Unable to determine the endianness of the binary file.")
