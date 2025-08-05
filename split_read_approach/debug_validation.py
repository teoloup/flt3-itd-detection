# Debug ITD validation logic
import logging

# Set up logging to see what's happening
logging.basicConfig(level=logging.DEBUG)

# Simulate the validation scenario from your log:
# - ITD candidate: 50bp at position 24 (original reference)
# - Reference creation: ITD 50bp inserted at position 76
# - Expected ITD region in modified reference: 76-126

print("=== ITD Validation Debug ===")
print("Scenario from your log:")
print("- ITD candidate: 50bp at position 24 (original reference)")
print("- ITD inserted at position 76 in modified reference")
print("- ITD region in modified reference: 76-126")
print()

# Simulate read alignment scenarios
print("Expected read behaviors:")
print()

print("1. WT reads aligned to ITD reference:")
print("   - Should align normally to positions 0-75")
print("   - Should have 50bp DELETION at positions 76-125 (missing ITD)")
print("   - Should continue aligning from position 126 onwards")
print("   - CIGAR example: 76M50D210M (76 matches, 50bp deletion, 210 matches)")
print()

print("2. ITD reads aligned to ITD reference:")
print("   - Should align normally to positions 0-75")  
print("   - Should MATCH perfectly at positions 76-125 (ITD region)")
print("   - Should continue aligning from position 126 onwards")
print("   - CIGAR example: 336M (all matches, no indels)")
print()

print("3. Current validation logic:")
print("   - Requires read to SPAN positions 76-126")
print("   - Counts MATCHES in ITD region")
print("   - If good matches → ITD supporting read")
print("   - If large deletion → WT read")
print()

print("4. Potential issues:")
print("   a) Read length: Are reads long enough to span the ITD region?")
print("   b) Alignment quality: Do reads actually align across the ITD?")
print("   c) Coordinate calculation: Are we looking at the right region?")
print()

# Simulate coordinate check
original_itd_pos = 24
insertion_pos = 76  
itd_length = 50
itd_start = insertion_pos
itd_end = insertion_pos + itd_length

print(f"Coordinate check:")
print(f"- Original ITD position: {original_itd_pos}")
print(f"- Insertion position in modified ref: {insertion_pos}")
print(f"- ITD region: {itd_start}-{itd_end}")
print(f"- Region length: {itd_end - itd_start}")
