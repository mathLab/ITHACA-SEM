<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady FRDG inviscidBurger MODIFIED, P=10</description>
    <executable>ADRSolver</executable>
    <parameters>InviscidBurger1D_FRDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">InviscidBurger1D_FRDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.28218</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.49972</value>
        </metric>
    </metrics>
</test>
