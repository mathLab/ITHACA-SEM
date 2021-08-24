<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady LFRHU explicit diffusion, regular quads, P=10, GAUSS Neumann bcs </description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_2D_LFRHU_regular_Neumann_GAUSS.xml</parameters>
    <files>
        <file description="Session File"> ExDiffusion_2D_LFRHU_regular_Neumann_GAUSS.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">6.83902e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.17219e-07</value>
        </metric>
    </metrics>
</test>




