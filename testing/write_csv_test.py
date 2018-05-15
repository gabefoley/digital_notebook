import csv
names= {"dog" : "harold", "cat" : "charlie"}

with open("./test_csv.csv", "w+") as file:
    for animal, name in names.items():
        file.write(name + ",")
        file.write("is a " + ",")
        file.write(animal)
