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
            set
            {
                Set(ref _selectedLayerParameters, value);
                if (value != null)
                {
                    SimulationL = value.L;
                    SimulationW = value.W;
                    SimulationH = value.H;
                    SelectedAcceleration = Accelerations[(int)value.AccelerationPattern];
                    SelectedCollision = Collisions[(int)value.CollisionMethod];
                }
                else
                {
                    SelectedAcceleration = "";
                    SelectedCollision = "";
                }
            }
        }

        public void AddLayerParameters()
        {
            LayerParameters layer = new LayerParameters();
            if (SimulationL != 0)
            {
                layer.L = SimulationL;
            }
            else
            {
                SimulationL = layer.L;
            }
            if (SimulationW != 0)
            {
                layer.W = SimulationW;
            }
            else
            {
                SimulationW = layer.W;
            }
            if (SimulationH != 0)
            {
                layer.H = SimulationH;
            }
            else
            {
                SimulationH = layer.H;
            }
            LayerParameters.Add(layer);
            SelectedLayerParameters = LayerParameters[LayerParameters.Count - 1];

            Debug.WriteLine("Parameter set added.");
        }

        public void RemoveLayerParameters()
        {
            if (SelectedLayerParameters != null)
            {
                LayerParameters.Remove(SelectedLayerParameters);
            }
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

        private ushort _simulationL;
        public ushort SimulationL
        {
            get => _simulationL;
            set
            {
                _simulationL = value;
                NotifyOfPropertyChange(() => SimulationL);
                if (LayerParameters.Count > 0)
                {
                    foreach (LayerParameters layer in LayerParameters)
                    {
                        layer.L = value;
                    }
                }

                // Remove once L and W are decoupled
                SimulationW = value;

                updateMemoryUsage();
            }
        }

        private ushort _simulationW;
        public ushort SimulationW
        {
            get => _simulationW;
            set
            {
                _simulationW = value;
                NotifyOfPropertyChange(() => SimulationW);
                if (LayerParameters.Count > 0)
                {
                    foreach (LayerParameters layer in LayerParameters)
                    {
                        layer.W = value;
                    }
                }
                updateMemoryUsage();
            }
        }

        private ushort _simulationH;
        public ushort SimulationH
        {
            get => _simulationH;
            set
            {
                _simulationH = value;
                NotifyOfPropertyChange(() => SimulationH);
                if (LayerParameters.Count > 0)
                {
                    foreach (LayerParameters layer in LayerParameters)
                    {
                        layer.H = value;
                    }
                }
                updateMemoryUsage();
            }
        }

        // Accleration control
        private BindableCollection<string> _accelerations = new BindableCollection<string>(new string[] { "None", "Accelerate", "Decelerate", "Bicone", "Hourglass" });
        public BindableCollection<string> Accelerations
        {
            get => _accelerations;
            set => Set(ref _accelerations, value);
        }

        private string _selectedAcceleration;
        public string SelectedAcceleration
        {
            get { return _selectedAcceleration; }
            set
            {
                _selectedAcceleration = value;
                NotifyOfPropertyChange(() => SelectedAcceleration);
                if (SelectedLayerParameters != null)
                {
                    if (value == Accelerations[0])
                    {
                        SelectedLayerParameters.AccelerationPattern = (Acceleration)0;
                    }
                    else if (value == Accelerations[1])
                    {
                        SelectedLayerParameters.AccelerationPattern = (Acceleration)1;
                    }
                    else if (value == Accelerations[2])
                    {
                        SelectedLayerParameters.AccelerationPattern = (Acceleration)2;
                    }
                    else if (value == Accelerations[3])
                    {
                        SelectedLayerParameters.AccelerationPattern = (Acceleration)3;
                    }
                    else if (value == Accelerations[4])
                    {
                        SelectedLayerParameters.AccelerationPattern = (Acceleration)4;
                    }
                }
            }
        }

        private BindableCollection<string> _collisions = new BindableCollection<string>(new string[] { "NN0", "NN1", "NN2", "NN3" });
        public BindableCollection<string> Collisions
        {
            get => _collisions;
            set => Set(ref _accelerations, value);
        }

        private string _selectedCollision;
        public string SelectedCollision
        {
            get { return _selectedCollision; }
            set
            {
                _selectedCollision = value;
                NotifyOfPropertyChange(() => SelectedCollision);
                if (SelectedLayerParameters != null)
                {
                    if (value == Collisions[0])
                    {
                        SelectedLayerParameters.CollisionMethod = (Collision)0;
                    }
                    else if (value == Collisions[1])
                    {
                        SelectedLayerParameters.CollisionMethod = (Collision)1;
                    }
                    else if (value == Collisions[2])
                    {
                        SelectedLayerParameters.CollisionMethod = (Collision)2;
                    }
                    else if (value == Collisions[3])
                    {
                        SelectedLayerParameters.CollisionMethod = (Collision)3;
                    }
                    updateMemoryUsage();
                }
            }
        }

        private System.UInt64 _estimatedMemoryUsage;
        public System.UInt64 EstimatedMemoryUsage
        {
            get => _estimatedMemoryUsage;
            set
            {
                _estimatedMemoryUsage = value;
                NotifyOfPropertyChange(() => EstimatedMemoryUsage);
            }
        }

        private void updateMemoryUsage()
        {
            int maximum_number_of_species = 0;
            uint maximum_reps = 0;
            bool collision_grid = false;
            foreach (LayerParameters layer in LayerParameters)
            {
                if (layer.Species.Count > maximum_number_of_species)
                {
                    maximum_number_of_species = layer.Species.Count;
                }

                if (layer.Repetitions > maximum_reps)
                {
                    maximum_reps = layer.Repetitions;
                }

                if (layer.CollisionMethod != (Collision)0)
                {
                    collision_grid = true;
                }
            }

            System.UInt64 grid_size = (System.UInt64)SimulationL * (System.UInt64)SimulationW * (System.UInt64)SimulationH; // simulation volume
            System.UInt64 memory = grid_size; // grid
            memory += grid_size * 4; // ordered
            if (collision_grid)
            {
                memory += grid_size; // collision grid
            }
            memory += 4 * maximum_reps * 3; // diffusion lengths
            memory += 2 * 27 * 3; // Surface3D.neighbors
            memory += grid_size; // Surface3D.adjacency
            memory += 8 * grid_size * (uint)maximum_number_of_species; // Surface3D.energy_grid
            memory += 4 * 27; // Surface3D.energies

            EstimatedMemoryUsage = memory >> 20;
            
        }
    }
}
