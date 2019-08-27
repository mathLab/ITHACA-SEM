<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Manufactured Compressible Poiseuille's flow to test IP</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>MMS_Compressible_Poiseuille_testIP.xml</parameters>
    <files>
        <file description="Session File">MMS_Compressible_Poiseuille_testIP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.71095e-06</value>
            <value variable="rhou" tolerance="1e-12">1.25767e-04</value>
            <value variable="rhov" tolerance="1e-12">3.98104e-04</value>
            <value variable="E" tolerance="1e-12">5.19340e-01</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">4.63915e-06</value>
            <value variable="rhou" tolerance="1e-12">5.07214e-04</value>
            <value variable="rhov" tolerance="1e-12">1.43833e-03</value>
            <value variable="E" tolerance="1e-12">1.26429e+00</value>
        </metric>
    </metrics>
</test>


