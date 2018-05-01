

import re
first = "AY753306.1:distance"
second = "AY753306.2:distance"

# print (first.split(".:")[1])
# print (second.split(".1:")[1])

print (re.split("[.]\d:", first)[1])
print (re.split("[.]\d:", second)[1])