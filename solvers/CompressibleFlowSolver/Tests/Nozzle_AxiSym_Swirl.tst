<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, axi-symmetric nozzle with swirl</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Nozzle_AxiSym_Swirl.xml</parameters>
    <files>
        <file description="Session File">Nozzle_AxiSym_Swirl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.05353</value>
            <value variable="rhou" tolerance="1e-12">2.33629</value>
            <value variable="rhov" tolerance="1e-12">102.209</value>
            <value variable="rhow" tolerance="1e-12">11.9746</value>
            <value variable="E" tolerance="1e-12">615434</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.25951</value>
            <value variable="rhou" tolerance="1e-12">2.79766</value>
            <value variable="rhov" tolerance="1e-12">60.5147</value>
            <value variable="rhow" tolerance="1e-12">34.8111</value>
            <value variable="E" tolerance="1e-12">260581</value>
        </metric>
    </metrics>
</test>


