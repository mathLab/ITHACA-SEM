<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Steady Linearised NavierStokes, 3D Soln with coupled solver, P=6</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>SinCos_LinNS_3DHom1D.xml</parameters>
    <files>
        <file description="Session File">SinCos_LinNS_3DHom1D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.42174e-15</value>
            <value variable="v" tolerance="1e-12">1.95785e-15</value>
            <value variable="w" tolerance="1e-12">8.16448e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.71606e-14</value>
            <value variable="v" tolerance="1e-12">1.77636e-14</value>
            <value variable="w" tolerance="1e-12">1.11177e-05</value>
        </metric>
    </metrics>
</test>
