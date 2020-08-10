<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Advection WDG MVM with --npz 3 </description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvection_WDG_3DHomo1D_MVM.xml --npz 3 </parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">UnsteadyAdvection_WDG_3DHomo1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.37343e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.50328e-05</value>
        </metric>
    </metrics>
</test>
