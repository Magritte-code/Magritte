// Component: CMake.Reference.DoxygenHTML

function Component()
{
    // Default constructor
}

Component.prototype.createOperations = function()
{
    // Create shortcut
    if (installer.value("os") === "win") {

        component.addOperation("CreateShortcut",
                               "@TargetDir@/doc/cmake-3.18/developer-reference/html/index.html",
                               "@StartMenuDir@/CMake Developer Reference.lnk");

    }

    // Call default implementation
    component.createOperations();
}
