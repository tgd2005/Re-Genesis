# 07/03/26

## Intuition behind weights

### wfile-1.txt

Exponential decrease with ~2:1 ratio of (0, 1) to (0, 0). For other states mindistance has no effect

### wfile-2.txt

Same as wfile-1 but same ratio of (0, 0) to (0, 0) - testing if that separation has any effect. (0, 0) = (0, 1) = 40562 which is pretty low.

### wfile-3.txt

Same as wfile-2 but (0, 0) = (0, 1) = 81623

### wfile-4.txt

Try bigger ratios between (0,0) and (0, 1) - now ~1:4.

### wfile-5.txt

Flip the ratio - (0, 1) is lower energy than (0, 0), so (0, 0):(0, 1) = 4:1
