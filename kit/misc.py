
def p2f(atm2bas_p):
    atm2bas_f = []
    for item in atm2bas_p:
        item_f = []
        for j in item:
            item_f.append(j+1)
        atm2bas_f.append(item_f)
    return atm2bas_f