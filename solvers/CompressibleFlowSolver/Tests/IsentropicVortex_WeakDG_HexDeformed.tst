<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex Hex Deformed WeakDG</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex_WeakDG_HexDeformed.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex_WeakDG_HexDeformed.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.00030573</value>
            <value variable="rhou" tolerance="1e-12">0.000809713</value>
            <value variable="rhov" tolerance="1e-12">0.000934604</value>
            <value variable="rhow" tolerance="1e-12">5.11694e-16</value>
            <value variable="E" tolerance="1e-12">0.00233073</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0002979</value>
            <value variable="rhou" tolerance="1e-12">0.00113882</value>
            <value variable="rhov" tolerance="1e-12">0.00124035</value>
            
            <value variable="rhow" tolerance="1e-12">7.87957e-15</value>
            <value variable="E" tolerance="1e-12">0.00271154</value>
        </metric>
    </metrics>
</test>


