# Orienteering Problem Instance Format

Each instance is stored as a plain text file with the following structure:

1. **First line:** Number of locations `n`.
2. **Second line:** Time limit (maximum allowed travel time).
3. **Next `n` lines:** Cost (distance) matrix of size `n × n`.  
   - Each line contains `n` integers separated by spaces.  
   - The diagonal entries are `0`.  
   - The matrix is typically symmetric, but asymmetric costs can also be used.
4. **Next `n` lines:** Scores for each location.  
   - One score per line.  
   - The order of scores corresponds to the order of locations in the matrix.
