using Caliburn.Micro;
using System.Diagnostics;
using BDDLL;

namespace BDSS.ViewModels
{
    public class MainViewModel : Screen
    {
        private BindableCollection<LayerParameters> _layerParameters = new BindableCollection<LayerParameters>();
        public BindableCollection<LayerParameters> LayerParameters
        {
            get => _layerParameters;
            set => Set(ref _layerParameters, value);
        }

        private LayerParameters _selectedLayerParameters;
        public LayerParameters SelectedLayerParameters
        {
            get => _selectedLayerParameters;
            set => Set(ref _selectedLayerParameters, value);
        }

        public void AddLayerParameters()
        {
            LayerParameters.Add(new LayerParameters());
            Debug.WriteLine("Parameter set added.");
        }

        public void RunAllFunctions()
        {
            if(LayerParameters.Count > 0)
            {
                Debug.WriteLine("Running simulation.");
                BallisticSimulation simulation = new BallisticSimulation();
                int ret = simulation.initializeSimulation(LayerParameters[0].L, LayerParameters[0].W, LayerParameters[0].H);
                if (ret != 0)
                {
                    return;
                }
                for (int i = LayerParameters.Count - 1; i >= 0; i--)
                {
                    var layerParameters = LayerParameters[i];
                    simulation.simulateFilm(layerParameters.Theta, layerParameters.L, layerParameters.H, layerParameters.Repetitions, layerParameters.Phi,
                        layerParameters.Turns, layerParameters.Seed, layerParameters.DiffusionSteps, layerParameters.Species, layerParameters.Weights,
                        layerParameters.AngularDeviation, layerParameters.PhiSweepPeriods, layerParameters.PhiDegree, layerParameters.ThetaSweep,
                        layerParameters.ThetaEnd, layerParameters.StepperResolution, layerParameters.MaterialSystem, layerParameters.AccelerationPattern,
                        layerParameters.CollisionMethod, layerParameters.DestinationDirectory);
                }
                simulation.cleanup();
            }
        }
    }
}
