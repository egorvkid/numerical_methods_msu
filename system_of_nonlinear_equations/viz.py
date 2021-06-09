import numpy as np
import matplotlib.pyplot as plt
x = np.arange(-1.0, 1.0, 0.01)
y = np.arange(0.0, 1.0, 0.01)
y1 = np.square(x)
y2 = np.sqrt(abs(1-np.square(x)))
plt.plot(x, y1, 'r', label='f1')
plt.plot(x, y2, 'b', label='f2')
plt.xlabel('x')
plt.ylabel('y')
plt.title('y(x)')
plt.grid()
plt.show()

x = np.arange(-1.0, 1.0, 0.01)
y = np.arange(0.0, 1.0, 0.01)

xg, yg = np.meshgrid(x, y)

zg = (xg**2 - yg)**2 + (xg**2 + yg**2 - 1)**2

plt.contour(xg, yg, zg, levels = 100)
plt.title("Level lines of squared discrepancy sum")
plt.colorbar()
plt.show()

with open("y0.txt", "r") as f:
    y = f.read().split()
with open("x0.txt", "r") as f:
    x = f.read().split()
with open("iter_num.txt", "r") as f:
    kol = f.read().split()
with open("ans_x.txt", "r") as f:
    ans_x = f.read().split()
with open("ans_y.txt", "r") as f:
    ans_y = f.read().split()
y = [float(i) for i in y]
x = [float(i) for i in x]
ans_y = [float(i) if i != '-nan(ind)' else 0 for i in ans_y]
ans_x = [float(i) if i != '-nan(ind)' else 0 for i in ans_x]
kol = [float(i) for i in kol]
masy = []
masx = []
for xi in x:
    for yi in y:
        masy.append(yi)
        masx.append(xi)
fig = plt.figure(figsize= (10,5))
plt.scatter(ans_x, ans_y)
plt.title("Simple_it")
plt.xlabel("x")
plt.ylabel("y")
plt.show()
plt.hexbin(masx, masy, kol, cmap = 'Accent')
plt.title("Simple it num iterations")
plt.xlabel("x_0")
plt.ylabel("y_0")
plt.colorbar()
plt.show()


with open("y01.txt", "r") as f:
    y = f.read().split()
with open("x01.txt", "r") as f:
    x = f.read().split()
with open("iter_num1.txt", "r") as f:
    kol = f.read().split()
with open("ans_x1.txt", "r") as f:
    ans_x = f.read().split()
with open("ans_y1.txt", "r") as f:
    ans_y = f.read().split()
y = [float(i) for i in y]
x = [float(i) for i in x]
ans_y = [float(i) if i != '-nan(ind)' else 0 for i in ans_y]
ans_x = [float(i) if i != '-nan(ind)' else 0 for i in ans_x]
kol = [float(i) for i in kol]
masy = []
masx = []
for xi in x:
    for yi in y:
        masy.append(yi)
        masx.append(xi)
fig = plt.figure(figsize=(10, 5))
plt.scatter(ans_x, ans_y)
plt.title("Newton")
plt.xlabel("x")
plt.ylabel("y")
plt.show()
plt.hexbin(masx, masy, kol, cmap = 'Accent')
plt.title("Newton num iterations")
plt.xlabel("x_0")
plt.ylabel("y_0")
plt.colorbar()
plt.show()

