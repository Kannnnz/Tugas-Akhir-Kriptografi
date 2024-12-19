import streamlit as st
import pandas as pd
import numpy as np

def import_sbox(file):
    try:
        df = pd.read_excel(file, header=None)
        return df.to_numpy().flatten()
    except Exception as e:
        st.error(f"Error: {e}")
        return None

def binary_representation(num, bits=8):
    return np.array([int(b) for b in format(num, f'0{bits}b')])

def optimized_walsh_hadamard(sbox, n=8, m=8):
    inputs = np.array([binary_representation(x, n) for x in range(2**n)])
    outputs = np.array([binary_representation(sbox[x], m) for x in range(2**n)])
    max_walsh = 0
    for u in range(1, 2**n):
        u_bin = binary_representation(u, n)
        for v in range(1, 2**m):
            v_bin = binary_representation(v, m)
            dot_u_x = (inputs @ u_bin) % 2
            dot_v_Sx = (outputs @ v_bin) % 2
            dot_result = (dot_u_x ^ dot_v_Sx)
            walsh_sum = np.sum(1 - 2 * dot_result)
            max_walsh = max(max_walsh, abs(walsh_sum))
    return 2**(n-1) - max_walsh / 2

def calculate_nonlinearity(sbox):
    return optimized_walsh_hadamard(sbox, n=8, m=8)

def calculate_sac_matrix(sbox):
    n = 8
    sac_matrix = np.zeros((n, n))
    for bit_input in range(n):
        for bit_output in range(n - 1, -1, -1):
            total_flips = 0
            for x in range(256):
                flipped_x = x ^ (1 << bit_input)
                original_output = (sbox[x] >> bit_output) & 1
                flipped_output = (sbox[flipped_x] >> bit_output) & 1
                total_flips += original_output ^ flipped_output
            sac_matrix[bit_input, n - 1 - bit_output] = total_flips / 256
    return np.round(sac_matrix, 5)

def calculate_sac_average(sac_matrix):
    return round(np.mean(sac_matrix), 5)

def calculate_bic_nl(sbox):
    n = 8
    total_nl = 0
    total_pairs = 0
    for bit1 in range(n):
        for bit2 in range(n):
            if bit1 != bit2:
                combined_sbox = [(sbox[x] >> bit1 & 1) ^ (sbox[x] >> bit2 & 1) for x in range(256)]
                nl = optimized_walsh_hadamard(combined_sbox, n=8, m=1)
                total_nl += nl
                total_pairs += 1
    return round(total_nl / total_pairs, 5)

def hamming_weight(x):
    return bin(x).count('1')

def calculate_bic_sac_matrix(sbox):
    n = len(sbox)
    bit_length = 8
    bic_sac_matrix = np.zeros((bit_length, bit_length))
    for i in range(bit_length):
        for j in range(i + 1, bit_length):
            independence_sum = 0
            for x in range(n):
                for bit_to_flip in range(bit_length):
                    flipped_x = x ^ (1 << bit_to_flip)
                    y1 = sbox[x]
                    y2 = sbox[flipped_x]
                    independence_sum += ((y1 >> i) & 1 ^ (y2 >> i) & 1) ^ ((y1 >> j) & 1 ^ (y2 >> j) & 1)
            bic_sac_matrix[i, j] = independence_sum / (n * bit_length)
            bic_sac_matrix[j, i] = bic_sac_matrix[i, j]
    return np.round(bic_sac_matrix, 6)

def calculate_bic_sac_average(bic_sac_matrix):
    bit_length = len(bic_sac_matrix)
    total_sum = np.sum(bic_sac_matrix) - np.trace(bic_sac_matrix)
    total_elements = bit_length * (bit_length - 1)
    return round(total_sum / total_elements, 6)


def calculate_lap(sbox):
    n = len(sbox)
    max_prob = 0
    for input_mask in range(1, 256):
        for output_mask in range(1, 256):
            count = 0
            for x in range(n):
                input_parity = bin(x & input_mask).count('1') % 2
                output_parity = bin(sbox[x] & output_mask).count('1') % 2
                if input_parity == output_parity:
                    count += 1
            prob = abs(count - n / 2) / n
            max_prob = max(max_prob, prob)
    return round(max_prob, 5)

def calculate_dap(sbox):
    n = len(sbox)
    max_prob = 0
    for input_diff in range(1, 256):
        for output_diff in range(256):
            count = 0
            for x in range(n):
                if x ^ input_diff < n and sbox[x] ^ sbox[x ^ input_diff] == output_diff:
                    count += 1
            prob = count / n
            max_prob = max(max_prob, prob)
    return round(max_prob, 5)

def save_to_excel(results, filename="hasil_kriptografi.xlsx"):
    with pd.ExcelWriter(filename) as writer:
        for name, matrix in results.items():
            if matrix is not None and len(matrix) > 0:
                pd.DataFrame(matrix).to_excel(writer, sheet_name=name, index=False, header=False)
    return filename

st.title("S-Box Cryptographic Strength Tester")

uploaded_file = st.file_uploader("Upload File S-Box (Excel)", type=["xlsx"])

if uploaded_file:
    sbox = import_sbox(uploaded_file)
    if sbox is not None:
        st.write("### S-Box yang Diimpor:")
        st.write(pd.DataFrame(sbox.reshape(16, 16)))

        operations = st.multiselect("Pilih Operasi Pengujian:",
                                    ["Nonlinearity (NL)", "Strict Avalanche Criterion (SAC)", "Bit Independence Criterion--Nonlinearity (BIC-NL)", "Bit Independence Criterion--Strict Avalanche Criterion (BIC-SAC)", "Linear Approximation Probability (LAP)", "Differential Approximation Probability (DAP)"])

        if st.button("Jalankan Pengujian"):
            results = {}
            for operation in operations:
                if operation == "Nonlinearity (NL)":
                    results["NL"] = [[calculate_nonlinearity(sbox)]]
                elif operation == "Strict Avalanche Criterion (SAC)":
                    sac_matrix = calculate_sac_matrix(sbox)
                    results["SAC Matrix"] = sac_matrix
                    sac_avg = calculate_sac_average(sac_matrix)
                    st.write(f"### SAC: {sac_avg}")
                elif operation == "Bit Independence Criterion--Nonlinearity (BIC-NL)":
                    results["BIC-NL"] = [[calculate_bic_nl(sbox)]]
                elif operation == "Bit Independence Criterion--Strict Avalanche Criterion (BIC-SAC)":
                    bic_sac_matrix = calculate_bic_sac_matrix(sbox)
                    results["BIC-SAC Matrix"] = bic_sac_matrix
                    bic_sac_avg = calculate_bic_sac_average(bic_sac_matrix)
                    st.write(f"### BIC-SAC: {bic_sac_avg}")
                elif operation == "Linear Approximation Probability (LAP)":
                    results["LAP"] = [[calculate_lap(sbox)]]
                elif operation == "Differential Approximation Probability (DAP)":
                    results["DAP"] = [[calculate_dap(sbox)]]

            for name, result in results.items():
                st.write(f"### {name}:")
                st.write(pd.DataFrame(result))

            file_path = save_to_excel(results)
            with open(file_path, "rb") as file:
                st.download_button("Download Hasil (.xlsx)", file, file_name="hasil_kriptografi.xlsx")
