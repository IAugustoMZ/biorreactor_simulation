import sys
import warnings
import matplotlib.pyplot as plt
sys.path.append('./')
from modules.biorreactor import BioReactor

# ignore warnings
warnings.filterwarnings('ignore')

# multiple constant
Cc0 = 1
Cs0 = 250.0
tb = 14

# class instance
br = BioReactor(Cc0=Cc0, Cs0=Cs0, production=1000.0, batch_time=tb)

# testing the Runge-Kutta's implementation
results = br.runge_kutta()

# create plot to show
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,3,1)
ax1.plot(results['batch_time'], results['cell_conc'], 'k-')
ax1.set_xlabel('Tempo de Batelada [h]', size = 14)
ax1.set_ylabel(r'Concentração de Células [$g.dm^{-3}$]', size = 14)
ax1.set_title('Crescimento Celular', size = 18)
ax1.grid(True)

ax2 = fig.add_subplot(1,3,2)
ax2.plot(results['batch_time'], results['substract_conc'], 'r-')
ax2.set_xlabel('Tempo de Batelada [h]', size = 14)
ax2.set_ylabel(r'Concentração de Substrato [$g.dm^{-3}$]', size = 14)
ax2.set_title('Consumo de Glicose', size = 18)
ax2.grid(True)

ax3 = fig.add_subplot(1,3,3)
ax3.plot(results['batch_time'], results['product_conc'], 'g-')
ax3.set_xlabel('Tempo de Batelada [h]', size = 14)
ax3.set_ylabel(r'Concentração de Produtos [$g.dm^{-3}$]', size = 14)
ax3.set_title('Produção de Etanol', size = 18)
ax3.grid(True)

br.calculate_volume()
print(br.volume)

plt.subplots_adjust(wspace=0.3)
plt.show()