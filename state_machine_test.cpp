#include <windows.h>

#include <iostream>
#include <stack>
using namespace std;

template <class EntityType, class MachineType, class Real=double>
class FSMSStateUpdateArgs
{
public:
	FSMSStateUpdateArgs(EntityType& entity, MachineType& machine, const Real& dt)
		: entity(entity), machine(machine), dt(dt)
	{}
	EntityType& entity;
	MachineType& machine;
	const Real& dt;
};

/**
 * @brief 有限状態マシン(Finite State Machine Static version)の状態クラス
 *
 */
template <class EntityType, class MachineType, class Real=double>
class FSMSState
{
public:
	typedef FSMSStateUpdateArgs<EntityType, MachineType, Real> FSMSStateUpdateArgsType;

	virtual ~FSMSState() = default;
	virtual void OnEnter(EntityType&) abstract;
	virtual void OnUpdate(FSMSStateUpdateArgsType&) abstract;
	virtual void OnExit(EntityType&) abstract;
};

/**
 * @brief 有限状態マシン(Finite State Machine Static version)の状態の空クラス
 *
 */
template <class EntityType, class MachineType, class Real=double>
class FSMSStateNull : public FSMSState<Real, EntityType, MachineType>
{
public:
	typedef FSMSStateUpdateArgs<EntityType, MachineType, Real> FSMSStateUpdateArgsType;

	virtual ~FSMSStateNull() = default;
	virtual void OnEnter(EntityType&) {}
	virtual void OnUpdate(FSMSStateUpdateArgsType&) {}
	virtual void OnExit(EntityType&) {}
};

/**
 * @brief 有限状態マシン(Finite State Machine Static version)の状態マシンクラス
 *
 */
template <class EntityType, class Real=double>
class FSMSMachine
{
public:
	typedef FSMSMachine<EntityType, Real> ThisType;
	typedef FSMSState<EntityType, ThisType, Real> StateType;
	typedef std::stack<StateType*> StateStackType;

public:
	FSMSMachine(ThisType* parent = nullptr) : _parent(parent)
	{} 

	void Update(EntityType& entity, Real dt)
	{
		if (_state_stack.empty()) return;

		StateType* state_top = _state_stack.top();

		typedef FSMSStateUpdateArgs<EntityType, ThisType, Real> FSMSStateUpdateArgsType;
		FSMSStateUpdateArgsType args(entity, *this, dt);
		state_top->OnUpdate(args);
	}

	void ChangeState(EntityType& entity, StateType& state)
	{
		while (!_state_stack.empty())
		{
			StateType* state_top = _state_stack.top();
			state_top->OnExit(entity);
			_state_stack.pop();
		}

		state.OnEnter(entity);
		_state_stack.push(&state);
	}

	void PushState(EntityType& entity, StateType& state)
	{
		state.OnEnter(entity);
		_state_stack.push(&state);
	}

	void PopState(EntityType& entity)
	{
		if (_state_stack.empty()) return;

		StateType* state_top = _state_stack.top();
		state_top->OnExit(entity);
		_state_stack.pop();
	}

	const StateType* CurrentState() const
	{
		if(_state_stack.empty()) return nullptr;

		return &_state_stack.top();
	}

	ThisType* Parent() const
	{
		return _parent;
	}

private:
	StateStackType _state_stack;
	ThisType* _parent;
};


// 以下、user code


class TestEntity;
class StateInit;
class StateStanding;
class StateShot;

class StateInitReady : public FSMSState<TestEntity, FSMSMachine<TestEntity> >
{
public:
	typedef StateInit ThisType;
	typedef TestEntity EntityType;
	typedef FSMSMachine<TestEntity> MachineType;

public:
	StateInitReady() {}
	virtual ~StateInitReady() {}
	virtual void OnEnter(EntityType&) override
	{
		cout << "StateInitReady::OnEnter" << endl;
	}

	virtual void OnUpdate(FSMSStateUpdateArgsType&) override;

	virtual void OnExit(EntityType&) override
	{
		cout << "StateInitReady::OnExit" << endl;
	}
};


class StateInit : public FSMSState<TestEntity, FSMSMachine<TestEntity> >
{
public:
	typedef StateInit ThisType;
	typedef TestEntity EntityType;
	typedef FSMSMachine<TestEntity> MachineType;

	MachineType machine;

	int count;

	StateInitReady mStateInitReady;

public:

	StateInit(EntityType& entity);
	virtual ~StateInit() {}
	virtual void OnEnter(EntityType&) override
	{
		cout << "StateInit::OnEnter" << endl;
	}
	
	virtual void OnUpdate(FSMSStateUpdateArgsType&) override;

	virtual void OnExit(EntityType&) override
	{
		cout << "StateInit::OnExit" << endl;
	}
};


class StateStanding : public FSMSState<TestEntity, FSMSMachine<TestEntity> >
{
public:
	typedef StateStanding ThisType;
	typedef TestEntity EntityType;
	typedef FSMSMachine<TestEntity> MachineType;

public:
	static ThisType& instance()
	{
		static ThisType a;
		return a;
	}

	virtual ~StateStanding() {}
	virtual void OnEnter(EntityType& entity) override;

	virtual void OnUpdate(FSMSStateUpdateArgsType&) override;

	virtual void OnExit(EntityType&) override
	{
		cout << "StateStanding::OnExit" << endl;
	}
};

class StateShot : public FSMSState<TestEntity, FSMSMachine<TestEntity> >
{
public:
	typedef StateShot ThisType;
	typedef TestEntity EntityType;
	typedef FSMSMachine<TestEntity> MachineType;

public:
	static ThisType& instance()
	{
		static ThisType a;
		return a;
	}

	virtual ~StateShot() {}
	virtual void OnEnter(EntityType&) override
	{
		cout << "StateShot::OnEnter" << endl;
	}

	virtual void OnUpdate(FSMSStateUpdateArgsType&) override;

	virtual void OnExit(EntityType&) override
	{
		cout << "StateShot::OnExit" << endl;
	}
};

class TestEntity
{
public:
	FSMSMachine<TestEntity> machine;

	bool is_shot;
	StateInit		mStateInit;
	StateStanding	mStateStanding;
	StateShot		mStateShot;


	TestEntity();


	void Update()
	{
		machine.Update(*this, 0.1);
	}
};

TestEntity::TestEntity() : mStateInit(*this)
{
	is_shot = false;
	machine.ChangeState(*this, this->mStateInit);
}


void StateInitReady::OnUpdate(FSMSStateUpdateArgsType& args)
{
	cout << "StateInitReady::OnUpdate " << endl;
	//args.machine.PopState(args.entity);
	auto Parent = args.machine.Parent();
	if (Parent)
	{
		args.machine.Parent()->ChangeState(args.entity, args.entity.mStateStanding);
	}
	//args.entity.machine.ChangeState(args.entity, args.entity.mStateStanding);
}


StateInit::StateInit(EntityType& entity) : count(0), machine(&entity.machine)
{
	machine.ChangeState(entity, this->mStateInitReady);
}

void StateInit::OnUpdate(FSMSStateUpdateArgsType& args)
{
	cout << "StateInit::OnUpdate " << ++count << endl;

	machine.Update(args.entity, args.dt);

	//args.machine.ChangeState(args.entity, args.entity.mStateStanding);
}

void StateStanding::OnEnter(EntityType& entity)
{
	cout << "StateStanding::OnEnter" << endl;
	entity.is_shot = false;
}

void StateStanding::OnUpdate(FSMSStateUpdateArgsType& args)
{
	cout << "StateStanding::OnUpdate" << endl;

	if (args.entity.is_shot)
	{
		args.machine.ChangeState(args.entity, args.entity.mStateInit);
	}
	else
	{
		args.machine.PushState(args.entity, args.entity.mStateShot);
	}
}

void StateShot::OnUpdate(FSMSStateUpdateArgsType& args)
{
	cout << "StateShot::OnUpdate" << endl;
	args.entity.is_shot = true;

	args.machine.PopState(args.entity);

	args.machine.Update(args.entity, args.dt);
}

int main()
{
	cout << "-- create a --" << endl;
	TestEntity a;

	cout << "-- create b --" << endl;
	TestEntity b;

	while(true)
	{
		cout << "-- a.Update() --" << endl;
		a.Update();
		cout << "-- b.Update() --" << endl;
		b.Update();
		Sleep(500);
	}

	return 0;
}