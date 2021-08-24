<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady LDG explicit diffusion, deformed quads, P=3 GAUSS</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_2D_LDG_deformed_GAUSS.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_2D_LDG_deformed_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00135534</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00191008</value>
        </metric>
    </metrics>
</test>




