<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Laminar Channel Flow 3DH1D, P=5, 8 Fourier modes (MVM) with flowrate control</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--npz 2 ChanFlow_3DH1D_Flowrate_MVM.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">ChanFlow_3DH1D_Flowrate_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">7.58667e-15</value>
            <value variable="p" tolerance="1e-8">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.76842e-18</value>
            <value variable="v" tolerance="1e-12">1.77358e-18</value>
            <value variable="w" tolerance="1e-12">2.44249e-14</value>
            <value variable="p" tolerance="1e-8">2.08802e-17</value>
        </metric>
    </metrics>
</test>
