for i in range(0, 17181, 500):
    final = i + 500 if (i + 500 < 17181)  else 17181
    print(i, final)

pct_dict = {"hops": 30, "hints": 35, "stropp": 21, "job": 56}

sorted_results = sorted(pct_dict.items(), key=lambda x: x[1], reverse=True)

for result in sorted_results:
    print (result)