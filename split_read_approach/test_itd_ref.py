# Test the new ITD reference creation logic
ref = 'AGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAAGGTACAG'

print(f'Original reference length: {len(ref)} bp')

# Simulate an ITD that duplicates region 100:150 (50bp)
dup_start = 100
dup_end = 150
duplicated_region = ref[dup_start:dup_end]

print(f'Duplicated region ({dup_start}:{dup_end}): {duplicated_region} ({len(duplicated_region)} bp)')

# Create ITD reference
itd_ref = ref[:dup_start] + duplicated_region + duplicated_region + ref[dup_end:]

print(f'ITD reference length: {len(itd_ref)} bp')
print(f'Length increase: {len(itd_ref) - len(ref)} bp')
print(f'Expected length increase: {len(duplicated_region)} bp')

# Show the structure
print(f'\nReference structure:')
print(f'Original: {ref[:20]}...{ref[-20:]}')
print(f'ITD:      {itd_ref[:20]}...{itd_ref[-20:]}')
print(f'\nDuplication check:')
print(f'Region appears twice: {duplicated_region in itd_ref}')
print(f'Region count in ITD ref: {itd_ref.count(duplicated_region)}')
