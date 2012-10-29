<test>
    <description>Channel Flow P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Test_ChanFlow_m3.xml</parameters>
    <files>
        <file description="Session File">Test_ChanFlow_m3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u">5.20394e-16</value>
            <value variable="v">0</value>
            <value variable="p">8.50369e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u">4.7462e-15</value>
            <value variable="v">3.04438e-16</value>
            <value variable="p">6.23945e-14</value>
        </metric>
    </metrics>
</test>


