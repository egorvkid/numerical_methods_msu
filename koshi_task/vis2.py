import matplotlib.pyplot as plt
def main():
	i = 0
	f = open('C:/Users/User/source/repos/CHW_8/CHW_8/y1.txt','r')
	l = []
	for line in f:
    		l.append(list(map(float, line.split())))
    		i+=1
	plt.plot(l[0],l[1],label = 'euler')
	plt.plot(l[0],l[2],label = 'euler_predict')
	plt.plot(l[0],l[3],label = 'runge2')
	plt.plot(l[0],l[4],label = 'runge4')
	plt.plot(l[0],l[5],label = 'adams3')
	plt.legend()
	plt.show()
main()