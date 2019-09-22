<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D Reaction Diffusion with bodyforce</description>
    <executable>ADRSolver</executable>
    <parameters>ReactionDiffusion2D.xml</parameters>
    <files>
        <file description="Session File">ReactionDiffusion2D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">1.41655e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">6.081e-07</value>
        </metric>
    </metrics>
</test>
