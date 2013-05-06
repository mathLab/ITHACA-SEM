<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D unsteady CG implicit diffusion, P = 5</description>
    <executable>ADRSolver</executable>
    <parameters>ImDiffusion_Hex_Periodic_m5.xml</parameters>
    <files>
        <file description="Session File">ImDiffusion_Hex_Periodic_m5.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00747775</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0052973</value>
        </metric>
    </metrics>
</test>
