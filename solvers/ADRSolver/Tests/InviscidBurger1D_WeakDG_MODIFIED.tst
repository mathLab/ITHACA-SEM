<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady WeakDG inviscidBurger MODIFIED, P=10</description>
    <executable>ADRSolver</executable>
    <parameters>InviscidBurger1D_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">InviscidBurger1D_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.28217</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.49986</value>
        </metric>
    </metrics>
</test>
