// Component: CMake.Documentation.SphinxHTML

function Component()
{
    // Default constructor
}

Component.prototype.createOperations = function()
{
    // Create shortcut
    if (installer.value("os") === "win") {

        component.addOperation("CreateShortcut",
                               "@TargetDir@/doc/cmake-3.18/html/index.html",
                               "@StartMenuDir@/CMake Documentation.lnk");

    }

    // Call default implementation
    component.createOperations();
}
