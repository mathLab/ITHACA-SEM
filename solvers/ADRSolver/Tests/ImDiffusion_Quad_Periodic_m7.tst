<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D unsteady CG implicit diffusion, P = 7</description>
    <executable>ADRSolver</executable>
    <parameters>ImDiffusion_Quad_Periodic_m7.xml</parameters>
    <files>
        <file description="Session File">ImDiffusion_Quad_Periodic_m7.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.53906e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-11">3.07765e-09</value>
        </metric>
    </metrics>
</test>
