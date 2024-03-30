using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
//using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

using BDDLL;

namespace BDSS
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void button_click(object sender, EventArgs e)
        {
            BDDLL.ManagedSimulationParameters parameters = new BDDLL.ManagedSimulationParameters();
            BDDLL.BallisticSimulation Wrapper = new BDDLL.BallisticSimulation();
            float theta = 85;
            List<byte> species = new List<byte> { 1 };
            List<List<float>> weights = new List<List<float>>{ };
            weights.Add(new List<float> { 1, 0 });
            weights.Add(new List<float> { 0, 1 });
            List<float> spread = new List<float> { 4, 2 };
            short[] inGrid = { 0 };
            short[] outGrid = { };
            Wrapper.simulateFilm(theta, 256, 256, 65536*4, 0, 0, 0, 5, species, weights, spread, 0, 0, false, theta, 0, "Si", 0, 0, "structures");
        }
    }
}
