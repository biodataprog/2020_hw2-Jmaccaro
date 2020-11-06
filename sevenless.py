#!/usr/bin/env python3

start = 0
end   = 99
divisor=7
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)

nums= []
for x in range(0,100):
        if (x%7 is not 0):
                nums.append(str(x))
print('\n'.join(nums))
