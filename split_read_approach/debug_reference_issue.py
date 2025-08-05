#!/usr/bin/env python3
"""
Debug script to analyze the reference construction issue
"""

# The reference from your file
reference_from_file = """AGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTAGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTATCTAGCTTTGCTGAATAATTAATCTCTGCTCATTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGAAAATATTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAAGGTACAG"""

# The ITD sequence from your VCF (316bp)
itd_sequence_from_vcf = """AGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTATCTAGCTTTGCTGAATAATTAATCTCTGCTCATTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGAAAATATTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAAT"""

# Expected insertion position (from VCF: position 160, but 0-based would be 159)
insertion_position = 160

print("=== DEBUG ANALYSIS ===")
print(f"Malformed reference length: {len(reference_from_file)}bp")
print(f"ITD sequence length: {len(itd_sequence_from_vcf)}bp")
print(f"Expected insertion position: {insertion_position}")

print("\n=== ISSUE DETECTION ===")

# Look for obvious duplications in the reference
print("Checking for duplicated sequences in the malformed reference...")

# Split into reasonable chunks to look for duplications
chunk_size = 50
seen_chunks = {}
duplicated_positions = []

for i in range(0, len(reference_from_file) - chunk_size + 1, 10):
    chunk = reference_from_file[i:i+chunk_size]
    if chunk in seen_chunks:
        duplicated_positions.append((i, seen_chunks[chunk], chunk))
        print(f"DUPLICATION found: positions {seen_chunks[chunk]} and {i}")
        print(f"  Sequence: {chunk[:30]}...")
    else:
        seen_chunks[chunk] = i

print(f"\nFound {len(duplicated_positions)} duplicated regions")

# Check if the ITD sequence is already present in the reference
print(f"\n=== ITD SEQUENCE ANALYSIS ===")
print(f"ITD sequence (first 50bp): {itd_sequence_from_vcf[:50]}")
print(f"ITD sequence (last 50bp): {itd_sequence_from_vcf[-50:]}")

# Check if ITD sequence is already present in reference
if itd_sequence_from_vcf in reference_from_file:
    pos = reference_from_file.find(itd_sequence_from_vcf)
    print(f"WARNING: ITD sequence already found in reference at position {pos}")
else:
    print("ITD sequence not found as complete match in reference")

# Check for partial matches
print("\nChecking for partial matches...")
for size in [100, 75, 50, 30]:
    itd_start = itd_sequence_from_vcf[:size]
    itd_end = itd_sequence_from_vcf[-size:]
    
    start_pos = reference_from_file.find(itd_start)
    end_pos = reference_from_file.find(itd_end)
    
    if start_pos >= 0:
        print(f"ITD start ({size}bp) found in reference at position {start_pos}")
    if end_pos >= 0:
        print(f"ITD end ({size}bp) found in reference at position {end_pos}")

print("\n=== WHAT SHOULD HAPPEN ===")
print("For a proper ITD insertion:")
print("1. Original reference should be clean FLT3 sequence")
print("2. ITD sequence should represent the duplicated portion")
print("3. Modified reference = original[0:pos] + ITD_sequence + original[pos:]")
print("4. No overlapping or duplicated content")

print("\n=== RECOMMENDATION ===")
print("The malformed reference suggests:")
print("1. The ITD sequence contains reference content that overlaps with existing sequence")
print("2. The insertion logic is creating duplicated content")
print("3. Need to fix either:")
print("   a) ITD sequence reconstruction (remove reference overlap)")
print("   b) Reference construction logic (handle overlaps)")

print("\n=== CRITICAL ANALYSIS: Do we lose legitimate ITDs? ===")

# Let's analyze what the ITD SHOULD represent vs what it actually represents
print("Understanding ITD biology:")
print("- ITD = Internal Tandem Duplication")
print("- A region of DNA gets duplicated and inserted")
print("- Example: ABCDEF → ABC[DEF]DEF (where DEF gets duplicated)")
print("")

# Let's try to reconstruct what the original reference should be
print("=== REVERSE ENGINEERING THE ORIGINAL REFERENCE ===")

# The malformed reference has the ITD sequence at position 159
# Let's see what the reference would look like without this
potential_original = reference_from_file[:159] + reference_from_file[159 + len(itd_sequence_from_vcf):]
print(f"Potential original reference length: {len(potential_original)}bp")
print(f"Original would be: {potential_original[:50]}...{potential_original[-50:]}")

# Now let's see if the ITD sequence makes sense
print(f"\n=== ITD SEQUENCE ANALYSIS ===")
print(f"ITD sequence length: {len(itd_sequence_from_vcf)}bp")

# Check if the ITD sequence contains parts that should be in the original reference
# vs parts that are truly duplicated/inserted content

# Compare ITD sequence to the potential original reference
print("\nChecking what parts of ITD sequence match the original reference:")

# Find the best alignment of ITD sequence to original reference
best_match_pos = -1
best_match_len = 0
best_overlap = 0

for start_pos in range(len(potential_original) - 20):
    for match_len in range(20, min(len(itd_sequence_from_vcf), len(potential_original) - start_pos) + 1):
        ref_segment = potential_original[start_pos:start_pos + match_len]
        if ref_segment in itd_sequence_from_vcf:
            if match_len > best_match_len:
                best_match_len = match_len
                best_match_pos = start_pos
                best_overlap = match_len

if best_match_pos >= 0:
    print(f"Best match: {best_match_len}bp of ITD matches original reference at position {best_match_pos}")
    matched_segment = potential_original[best_match_pos:best_match_pos + best_match_len]
    itd_match_start = itd_sequence_from_vcf.find(matched_segment)
    print(f"  Reference segment: {matched_segment[:30]}...")
    print(f"  Found in ITD at position: {itd_match_start}")
    
    # Check if there are novel insertions in the ITD
    itd_before_match = itd_sequence_from_vcf[:itd_match_start]
    itd_after_match = itd_sequence_from_vcf[itd_match_start + best_match_len:]
    
    print(f"\nITD components:")
    print(f"  Before match: {len(itd_before_match)}bp - {itd_before_match[:30]}...")
    print(f"  Matched region: {best_match_len}bp - (reference duplication)")
    print(f"  After match: {len(itd_after_match)}bp - {itd_after_match[:30]}...")
    
    # Check if before/after parts are also in reference (indicating more duplication)
    if itd_before_match and itd_before_match in potential_original:
        before_pos = potential_original.find(itd_before_match)
        print(f"  WARNING: 'Before match' part also found in reference at position {before_pos}")
    
    if itd_after_match and itd_after_match in potential_original:
        after_pos = potential_original.find(itd_after_match)
        print(f"  WARNING: 'After match' part also found in reference at position {after_pos}")

print(f"\n=== LEGITIMATE ITD DETECTION ===")
print("Types of legitimate ITDs we need to preserve:")
print("1. Pure duplications: existing sequence gets duplicated")
print("2. Duplications with novel insertions: existing sequence + new DNA")
print("3. Complex rearrangements: multiple duplicated segments + insertions")
print("")
print("Our fix should handle:")
print("✓ Case 1: Detect when ITD = reference duplication, avoid double-insertion")
print("? Case 2: Preserve novel insertions while avoiding reference duplication")
print("? Case 3: Handle complex ITDs correctly")

print(f"\n=== RISK ASSESSMENT ===")
risk_factors = []
if len(itd_sequence_from_vcf) > 200:
    risk_factors.append("Large ITD (>200bp) - likely contains reference + novel content")
if itd_sequence_from_vcf in reference_from_file:
    risk_factors.append("Complete ITD sequence found in reference - high duplication risk")

novel_content_estimate = 0
if best_match_pos >= 0:
    total_matched = best_match_len
    novel_content_estimate = len(itd_sequence_from_vcf) - total_matched
    if novel_content_estimate > 20:
        risk_factors.append(f"Estimated {novel_content_estimate}bp novel content - risk of losing insertions")

print(f"Risk factors identified: {len(risk_factors)}")
for risk in risk_factors:
    print(f"  - {risk}")

if novel_content_estimate > 20:
    print(f"\n*** WARNING: This ITD may contain {novel_content_estimate}bp of novel insertions ***")
    print("Current fix might trim legitimate novel content!")
    print("Need more sophisticated approach to separate:")
    print("  1. Reference duplications (should not be double-inserted)")
    print("  2. Novel insertions (should be preserved)")
else:
    print(f"\nLow risk: ITD appears to be mostly reference duplication")
    print("Current fix should work correctly")
