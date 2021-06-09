import matplotlib.pyplot as plt
def main():
	i = 0
	f = open('C:/Users/User/source/repos/CHW_9/CHW_9/nonlinear.txt','r')
	l = []
	for line in f:
    		l.append(list(map(float, line.split())))
    		i+=1
	plt.plot(l[0],l[1],label = 'newton')
	plt.plot(l[0],l[2],label = 'chordes')
	plt.legend()
	plt.show()
main()