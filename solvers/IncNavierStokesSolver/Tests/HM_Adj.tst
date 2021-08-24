<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fourier Half Mode Adjoint Basis, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>HM_Adj.xml</parameters>
    <files>
        <file description="Session File">HM_Adj.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.24579e-16</value>
            <value variable="v" tolerance="1e-12">2.7432e-16</value>
            <value variable="w" tolerance="1e-12">2.21876e-16</value>
            <value variable="p" tolerance="1e-12">1.20662e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.77636e-15</value>
            <value variable="v" tolerance="1e-12">1.08247e-15</value>
            <value variable="w" tolerance="1e-12">8.32667e-16</value>
            <value variable="p" tolerance="1e-12">4.69775e-14</value>
        </metric>
    </metrics>
</test>


