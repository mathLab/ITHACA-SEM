<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D homogeneous 1D Channel Flow, SEM parallelisation (2 proc)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_Parallel_mode1.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">ChanFlow_3DH1D_Parallel_mode1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">3.44241e-14</value>
            <value variable="v" tolerance="1e-6">2.09484e-13</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">4.8457e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">2.72504e-13</value>
            <value variable="v" tolerance="1e-6">3.29175e-13</value>
            <value variable="w" tolerance="1e-6">1.22076e-18</value>
            <value variable="p" tolerance="1e-6">1.33911e-10</value>
        </metric>
    </metrics>
</test>


