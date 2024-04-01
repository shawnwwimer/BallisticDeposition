using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Caliburn.Micro;
using System.Windows;
using System.Diagnostics;
using BDSS.ViewModels;
using BDSS.Views;

namespace BDSS
{
    public class AppBootstrapper : BootstrapperBase
    {
        public AppBootstrapper()
        {
            Initialize();
        }

        protected override async void OnStartup(object sender, StartupEventArgs e)
        {
            await DisplayRootViewForAsync<BDSS.ViewModels.MainViewModel>();
        }

        //protected override void Configure()
        //{
        //    ViewLocator.LocateForModelType += (modelType, displayLocation, context) =>
        //    {
        //        Debug.WriteLine($"Trying to locate view for {modelType} in location {displayLocation}.");
        //        return ViewLocator.LocateForModelType(modelType, displayLocation, context);
        //    };
        //}
    }
}
