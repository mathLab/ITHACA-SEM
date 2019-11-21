<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Laminar Channel Flow 3DH1D, P=4, 12 Fourier modes (MVM) with flowrate control and forcing direction explicitly defined</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_FlowrateExplicit_MVM.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_3DH1D_FlowrateExplicit_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">2.58781e-08</value>
            <value variable="v" tolerance="1e-8">1.0607e-15</value>
            <value variable="w" tolerance="1e-8">0</value>
            <value variable="p" tolerance="1e-8">5.32757e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">4.57503e-08</value>
            <value variable="v" tolerance="1e-8">7.23752e-15</value>
            <value variable="w" tolerance="1e-8">8.74374e-25</value>
            <value variable="p" tolerance="1e-8">7.74547e-13</value>
        </metric>
    </metrics>
</test>
