<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Advection FRSD MVM</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvection_FRSD_3DHomo1D_MVM.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvection_FRSD_3DHomo1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.65576e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6.23615e-05</value>
        </metric>
    </metrics>
</test>
