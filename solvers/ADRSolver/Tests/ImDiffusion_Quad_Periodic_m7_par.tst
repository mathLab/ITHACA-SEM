<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D unsteady CG implicit diffusion, P = 7</description>
    <executable>ADRSolver</executable>
    <parameters>--use-scotch ImDiffusion_Quad_Periodic_m7.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">ImDiffusion_Quad_Periodic_m7.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">1.71218e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">5.10955e-08</value>
        </metric>
    </metrics>
</test>

