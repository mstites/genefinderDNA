#For testing items for gene_finder

string = "hello world2 world1"
print(string)
test = string[1]
print(test)
print(type("A"))
string_n = string.replace("world","")
print(string_n)


def find(word, letter):
    index = 0
    while index < len(word):
        if word[index] == letter:
            return index
        index = index + 1
    return -1
word = 'banana'
count = 0
for letter in word:
    if letter == 'a':
        count = count + 1
print(count)

for i in range(len(numbers)):
    numbers[i] = numbers[i] * 2
