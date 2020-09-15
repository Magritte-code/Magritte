// Component: CMake

function Component()
{
    // Do not show component selection page
    installer.setDefaultPageVisible(QInstaller.ComponentSelection, false);
}

Component.prototype.createOperations = function()
{
    // Create shortcut
    if (installer.value("os") === "win") {



        component.addOperation("CreateShortcut",
                               "@TargetDir@/doc/cmake-3.18/cmake.org.html",
                               "@StartMenuDir@/CMake Web Site.lnk");

        component.addOperation("CreateShortcut",
                               "@TargetDir@/cmake-maintenance.exe",
                               "@StartMenuDir@/CMake Maintenance Tool.lnk");
    }

    // Call default implementation
    component.createOperations();
}
