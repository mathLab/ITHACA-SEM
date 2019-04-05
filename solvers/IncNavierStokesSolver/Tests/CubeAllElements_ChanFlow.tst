<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, all elements, HDF5 input, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CubeAllElements_ChanFlow.xml</parameters>
    <files>
        <file description="Session File">CubeAllElements_ChanFlow.xml</file>
        <file description="Geometry File">CubeAllElements_ChanFlow.nekg</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.3823e-15</value>
            <value variable="v" tolerance="1e-12">1.51882e-15</value>
            <value variable="w" tolerance="1e-12">6.24669e-15</value>
            <value variable="p" tolerance="1e-8">1.19474e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">7.06676e-15</value>
            <value variable="v" tolerance="1e-12">8.366e-15</value>
            <value variable="w" tolerance="1e-12">6.00076e-14</value>
            <value variable="p" tolerance="1e-8">1.49414e-12</value>
        </metric>
    </metrics>
</test>
