<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Advection FRHU MVM</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvection_FRHU_3DHomo1D_MVM.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvection_FRHU_3DHomo1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">6.34421e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">3.96748e-05</value>
        </metric>
    </metrics>
</test>
