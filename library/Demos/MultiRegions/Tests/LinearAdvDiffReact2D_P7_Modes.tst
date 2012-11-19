<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Steady ADR 2D problem with CG P=7</description>
    <executable>SteadyAdvectionDiffusionReaction2D</executable>
    <parameters>LinearAdvDiffReact2D_P7_Modes.xml</parameters>
    <files>
        <file description="Session File">LinearAdvDiffReact2D_P7_Modes.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0.00890901</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.0102694</value>
        </metric>
    </metrics>
</test>


