def quickSort(lst):
    if len(lst) <= 1: 
        return lst
    smaller = [x for x in lst[1:] if x < lst[0]]
    larger = [x for x in lst[1:] if x >= lst[0]]
    return quickSort(smaller) + [lst[0]] + quickSort(larger)

# Main Function
if __name__ == '__main__':    
    lst = [2,4,5,1]
    print quickSort(lst)    
    lst2 = [2,4,5,1,20,-1,999,-50,-25,25,25,25]
    print quickSort(lst2)    
