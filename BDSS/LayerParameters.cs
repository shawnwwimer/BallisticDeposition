using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Caliburn.Micro;

namespace BDSS
{
    public class LayerParameters : PropertyChangedBase
    {
        public LayerParameters()
        {
            Theta = 85;
            L = 256;
            W = L;
            H = 256;
            Repetitions = 65536 * 4;
            Phi = 0;
            Turns = 0;
            Seed = 0;
            DiffusionSteps = 5;
            Species = new List<byte> { 1 };
            Weights = new List<List<float>> { };
            Weights.Add(new List<float> { 1, 0 });
            Weights.Add(new List<float> { 0, 1 });
            AngularDeviation = new List<float> { 0f, 0f };
            PhiSweepPeriods = 0;
            PhiDegree = 0;
            ThetaSweep = false;
            ThetaEnd = Theta;
            StepperResolution = 0;
            MaterialSystem = "Si";
            AccelerationPattern = (Acceleration)1;
            CollisionMethod = (Collision)1;
            DestinationDirectory = "structures/";
        }

        private float _theta;
        public float Theta 
        { 
            get => _theta;
            set
            {
                _theta = value;
                NotifyOfPropertyChange(() => Theta);
            }
        }
        private UInt16 _l;
        public UInt16 L 
        { 
            get => _l;
            set
            {
                _l = value;
                W = value;
                NotifyOfPropertyChange(() => L);
            }
        }
        private UInt16 _w;
        public UInt16 W 
        { 
            get => _w;
            set
            {
                _w = value;
                NotifyOfPropertyChange(() => W);
            }
        }
        public UInt16 H { get; set; }
        public UInt32 Repetitions { get; set; }
        private float _phi;
        public float Phi 
        { 
            get => _phi;
            set
            {
                _phi = value;
                NotifyOfPropertyChange(() => Phi);
            }
        }
        public float Turns { get; set; }
        public uint Seed { get; set; }
        public UInt16 DiffusionSteps { get; set; }
        public List<byte> Species { get; set; }
        public List<List<float>> Weights { get; set; }
        public List<float> AngularDeviation { get; set; }
        public int PhiSweepPeriods { get; set; }
        public float PhiDegree { get; set; }
        public bool ThetaSweep { get; set; }
        public float ThetaEnd { get; set; }
        public UInt32 StepperResolution { get; set; }
        private String _material_system;
        public String MaterialSystem 
        { 
            get => _material_system;
            set
            {
                _material_system = value;
                NotifyOfPropertyChange(() => MaterialSystem);
            }
        }
        public Acceleration AccelerationPattern { get; set; }
        public Collision CollisionMethod { get; set; }
        public String DestinationDirectory { get; set; }
    }
}
