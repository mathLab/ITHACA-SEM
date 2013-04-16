<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady LFRHU explicit diffusion, regular quads, P=10, Neumann bcs</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_2D_LFRHU_regular_Neumann.xml</parameters>
    <files>
        <file description="Session File"> ExDiffusion_2D_LFRHU_regular_Neumann.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.88972e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">7.42037e-08</value>
        </metric>
    </metrics>
</test>




