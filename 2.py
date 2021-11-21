from random import randint
def main():
    p, j, D = CreateEC()
    G = Graph()
    n = 0
    # вычисление степени изогений
    for l in Primes():
        if kronecker_symbol(D, l) == 1:
            list_j = j_invariants(p, j, l)
            graph(list_j, G)
            n += 1
        if n == 5:
            break
    paintG = G.graphplot(vertex_labels = True, vertex_size = 1000, edge_thickness = 0.25, edge_color = 'blue', layout = 'circular', vertex_color = 'magenta', graph_border = True)
    paintG.show(figsize=[15,15])

def graph(j_inv, G_star):
    G = Graph()
    i = j_inv[len(j_inv) - 1]
    for j in j_inv:
        G.add_edge(i, j)
        G_star.add_edge(i, j)
        i = j
    paintG = G.graphplot(vertex_labels = True, vertex_size = 1000, edge_thickness = 0.5, edge_color = 'blue', layout = 'circular', vertex_color = 'yellow', graph_border = True)
    paintG.show(figsize=[15,15])

# вычисление перестановки j-инвариантов изогенных эллиптических кривых
def j_invariants(p, j, l):
    # генерация конечного поля характеристики p
    K = GF(p)
    # возвращает базу данных классических модульных функций, то есть многочленов Phi_li (X, Y)
    PHI = ClassicalModularPolynomialDatabase()
    # приведение базы полиномов в поле K
    phi = (PHI[l]).change_ring(K)
    # строит кольцо полиномов
    R.<x> = PolynomialRing(K)
    j_inv = [j]
    inv = j
    # решаем уравнение от j и переменной x корнями будут инварианты
    while True:
        # вычисляем 2 корня
        r = phi(inv, x).roots(ring = K, multiplicities = False)
        if j_inv.count(r[0]) != 0:
            if j_inv.count(r[1]) != 0:
                break
            else:
                # выбираем кривую с j = r[1]
                j_inv.append(r[1])
                inv = r[1]
        else:
            # выбираем кривую с j = r[0]
            j_inv.append(r[0])
            inv = r[0]
    print ('l = {}: {}'.format(l, j_inv))
    return j_inv

# вычисление параметров и генерация кривой
def CreateEC():
    p, cl, EC, T, D = 0, 0, 0, 0, 0
    while(1):
        # характеристика поля
        p = randint(10**4, 10**5)
        p = next_prime(p)
        # ненулевые коэффициенты кривой
        A, B = randint(1, p), randint(1, p)
        # генерация конечного поля характеристики p
        K = GF(p)
        # генерация кривой над конечным полем K
        E = EllipticCurve(K, [A, B])
        # вычисление значение следа T эндоморфизма Фробениуса
        T = E.trace_of_frobenius()
        # вычисление дискриминанта
        D = T**2-4*p
        # вычисление числа классов
        cl = QuadraticField(D).class_number()
        if(cl in Primes() and cl > 30 and cl < 100):
            break
    # вычисление полинома Гильберта для дискриминанта D
    H = hilbert_class_polynomial(D)
    # j инвариант кривой
    j = H.roots(ring = GF(p), multiplicities = False)
    if (j == []):
        print ("Error j")
        CreateEC()
    print ('p = {}, A = {}, B = {}: {}'.format(p, A, B, E))
    print ('D = {}, number of class = {}'.format(D, cl))
    print ('j-invariant = {}\n'.format(j[0]))
    return p, j[0], D

if __name__ == "__main__":
    main()