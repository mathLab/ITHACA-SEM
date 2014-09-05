<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady LFRHU explicit diffusion, deformed quads, P=3 GAUSS</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_2D_LFRHU_deformed_GAUSS.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_2D_LFRHU_deformed_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0012858</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0029901</value>
        </metric>
    </metrics>
</test>




